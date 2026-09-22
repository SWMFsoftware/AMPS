# UStandaloneProducts — Step 7 contract tests

Run from the repository root:

```bash
./srcEarth/test/UStandaloneProducts/run_test.sh
```

The dependency-free C++11 suite exercises the same `StandaloneProductContract.h` used
by gridless and standalone Mode3D dispatch. It verifies:

- aliases for DIPOLE, IGRF, T96, T01, T05/TS05, TA15N/B and TA16;
- cutoff-only, flux/spectrum-only and combined product selection;
- POINTS, TRAJECTORY and SHELLS domain validation;
- one authoritative epoch for field drivers, boundary spectrum, ephemeris and output;
- fail-fast behavior for unsupported models and unvalidated external-driver columns or
  units; and
- a machine-readable manifest containing model, representation, epoch and products.

This suite provides numerical/logic contract coverage without model libraries. The
configured I-F01–I-F07 and C1–C19 profiles remain responsible for linked Geopack model
values, driver-table interpolation, gridless/mesh agreement and observational data.
