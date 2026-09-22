# UBoundaryProducts — Step 6 reference tests

Run from the repository root:

```bash
./srcEarth/test/UBoundaryProducts/run_test.sh
```

The suite compiles the production header-only boundary/product kernel and production
`cSpectrum` evaluator with C++11. A second C++17 target compiles and links the production
`gridless/AnisotropicSpectrum.cpp` adapter, preventing source-level declaration/order
regressions from being hidden by helper-only tests. It requires no AMPS executable, MPI,
SPICE, Geopack, or SWMF libraries.

These are reference comparisons, not only “finite/nonzero” smoke checks:

| ID | Comparison reference |
|---|---|
| U-F01 | closed-form POWER_LAW, POWER_LAW_CUTOFF, LIS force-field and Band equations; exact log-log table interpolation; per-MeV/nucleon ↔ per-particle-J round trip |
| U-F05 | analytically clipped constant-integrand energy channel |
| U-F07 | beta-function means for `sin²(alpha)` and `cos²(alpha)`; normalized sphere/hemisphere mean exactly one |
| U-F07 adapter | production `EvalAnisotropyFactor`: raw and unit-mean `sin²`/`cos²`, normalized day/night factors, NaN fallback, and invalid-model rejection |
| U-F08 | geometric-mean log-intensity interpolation, exact-row selection, gap tie-breaking, ZERO/FAIL policies |
| U-F09 | analytic top-hat and triangular detector folds and isotropic `F_planar=F_omni/4` |
| F1/F16 kernel | fully allowed/blocked limiting behavior |
| F4 kernel | reintegration of the emitted differential-spectrum samples reproduces the reported integral |
| Step-6 phase space | direct `j_local=A(p_local²/p_boundary²)j_boundary` value |

The linked F1, F2, F4, F5, F11, F12, F15 and F16 cases remain the end-to-end tests of
trajectory access and output files in a configured AMPS build. This unit suite isolates
the boundary distribution, units, interpolation, uncertainty, and product quadrature so
a field or trajectory failure cannot masquerade as a flux-postprocessing defect.
