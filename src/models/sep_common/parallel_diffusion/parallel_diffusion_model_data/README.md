# Parallel-diffusion companion data, bundle revision 1.4

This directory is the numerical reference bundle for
`PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md`, revision 1.4 (Section 22). It is test and
provenance data for a future implementation, not the diffusion-coefficient library.

## Provenance

- `NLGCE_F_2014_parallel.csv` and `NLGCE_F_2014_perpendicular.csv` are the coefficients
  d_ijkl of Qin, G. & Zhang, L.-H. (2014), ApJ 787, 12, doi:10.1088/0004-637X/787/1/12,
  Tables 3 and 4 (arXiv:1401.1950v2). They are byte-identical to the files whose SHA-256
  digests were published in revisions 1.2 and 1.3 of the specification:

      7cdc5ab9cddda0c7295a25bfcaff4ba3422002da12ac1b4660a98dcc64eaa762  NLGCE_F_2014_parallel.csv
      7371b6557f40c5efa7a48460a5f1b2971374a2c43d0eb3e209f0a8e62670356a  NLGCE_F_2014_perpendicular.csv

  In revision 1.4 all 576 values were compared again, one by one, with the arXiv text of
  both tables. The perpendicular row (j,k,l)=(0,0,1) is garbled in the text extraction of
  the v2 PDF; it was confirmed from arXiv:1401.1950v1, where it is printed cleanly.
- All other files were regenerated in revision 1.4 by an independent implementation of the
  specification equations. The original revision 1.2/1.3 versions of these files could not
  be located. Every value printed in the specification (Sections 11.4 and 15) agrees with
  the regenerated values to its printed precision, and the 300-state audit sample was
  reproduced exactly from its stated generator and seed.

## Files

| File | Contents |
|---|---|
| `NLGCE_F_2014_parallel.csv`, `NLGCE_F_2014_perpendicular.csv` | 48 rows each: `j,k,l,d_i0,...,d_i5`. Build d[i,j,k,l] from the labels; do not infer a memory layout from the row order. |
| `benchmark_points.json` | Full-precision fixtures and their comparison tolerances. Particle kinematics: proton, electron, and alpha particle given as energy per nucleon. Spectrum and pitch-angle constants. Exact QLT, Eq. (30), with the inertial approximation, Eq. (31). NLGCE-F values and analytical log-derivatives at states A–G. NLGCE-N values, a_x, a'^2, and finite-difference log-derivatives. NLGC-E values. NLPA with a supplied perpendicular coefficient. Lorentzian and Gaussian broadened slab. |
| `fit_error_samples.csv` | All 300 audit states: the drawn natural-log inputs, the ratios, a_x, a'^2, both backend pairs, signed discrepancies, residuals, and refinement checks. |
| `fit_error_summary.json` | Generator, conventions, solver controls, environment, the six-bin statistics, overall quantiles and maxima, and the worst state. |
| `configuration_examples.json` | The non-normative templates of Section 14.6. Placeholders only; it contains no physical defaults. |
| `reference_verification.py` | Standalone verifier. It recomputes every fixture from the equations and checks the digests. |
| `SHA256SUMS` | Digests of every other file in this directory. |

## Commands

From this directory:

```sh
python3 reference_verification.py              # digests, schema, all deterministic fixtures (~1 min)
python3 reference_verification.py --audit      # also re-solve the 300-state audit and print the six bins
python3 reference_verification.py --broadened  # also recompute the broadened-slab fixtures (a few minutes)
```

Dependencies are Python 3.9+, NumPy and SciPy. `fit_error_summary.json` records the
generation environment. A different compatible NumPy/SciPy may change the last
floating-point digits, but not beyond the declared tolerances. The `--audit` draw relies on
NumPy's `default_rng` (PCG64) stream, which NumPy keeps stable across versions.

## Conventions

- Natural logarithms throughout. The NLGCE-F inputs are the original ratios r_L/ell_s,
  f_s = dB_s^2/(dB_s^2+dB_2^2), eps^2 = total variance/B0^2, and ell_s/ell_2, with ell_a the
  spectral bend-over lengths.
- a_x follows Eq. (45): its denominator is (xi/(1+xi))/eps + eps/(2 xi), with
  xi = (r_L/ell_s)/(C(nu) eps) and eps = sqrt(total variance)/B0. The expression
  `(xi/(1+xi))/epsilon + epsilon/(2*xi)` is that denominator, not a_x itself.
- Dimensionless fixtures use ell_s = 1, v = 1, B0 = 1. To drive an SI interface, pick any
  particle and B0, set ell_s = r_L/(r_L/ell_s), ell_2 = ell_s/(ell_s/ell_2),
  dB^2 = eps^2 B0^2, and compare lambda/ell_s.
- Every tolerance is a numerical comparison control. None is a physical accuracy statement.
  Agreement between the polynomial and the integral closure is a separate measured
  quantity (Section 11.4).
