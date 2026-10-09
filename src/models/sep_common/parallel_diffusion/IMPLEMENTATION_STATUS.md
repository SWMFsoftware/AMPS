# Parallel diffusion implementation status

Updated: 2026-10-08. Working tree: uncommitted; no commit or push performed.

## Input contract

- Scientific specification: `PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md`,
  revision 1.4.
- Companion bundle: `parallel_diffusion_model_data/`, bundle revision 1.4.
- Published NLGCE-F array identities:
  `7cdc5ab9cddda0c7295a25bfcaff4ba3422002da12ac1b4660a98dcc64eaa762`
  (parallel) and
  `7371b6557f40c5efa7a48460a5f1b2971374a2c43d0eb3e209f0a8e62670356a`
  (perpendicular).
- Approved scientific departures: none.
- User decision: D14 remains deferred. No srcSEP3D/srcSEP parser or mover
  binding is part of this change.

## Repository decisions

- Reusable API and implementation:
  `src/models/sep_common/parallel_diffusion/`.
- Language/dependencies: C++17 standard library only; fixture orchestration is
  Python 3 standard library. The independent bundle verifier additionally
  requires NumPy/SciPy.
- Standalone archive: `libparallel_diffusion.a`.
- NLGCE-F generation: `generate_nlgce_coefficients.py` validates exact CSV
  index coverage and emits `nlgce_f_coefficients.inc` without refitting or
  changing literals.
- The established seven-member `sep_common.a` is unchanged until PD11 updates
  both application ownership audits in one qualified host-integration change.

## Stage state

| Stage | State | Evidence or remaining boundary |
| --- | --- | --- |
| PD00 | complete | Bundle digests and default/audit/broadened verifier modes pass. |
| PD01 | complete | SI state, statuses, registry, provenance, fingerprints, transactional manager and active function pointer tested. |
| PD02 | complete | Five explicit prescriptions and their identities/limits tested. |
| PD03 | complete | Canonical supplied-spectrum representation, all named Section 8.7 sidedness/component/unit/frozen-flow conversions, named tail/coverage policies, normalized smooth and Equation (35) multirange spectra (including the `s=1` limit), break-partitioned adaptive integration and immutable-state contract implemented. |
| PD04 | complete | Pitch-angle conversion, prescribed shape, exact smooth-spectrum QLT, multirange/supplied-spectrum routes and inertial approximation implemented. Exact Eq. (30), independent multirange D-mu-mu, and divergent dissipation-tail checks pass. |
| PD05 | complete | Both exact published arrays compiled from audited CSV, fit box enforced, Horner values and analytic derivatives implemented; A–G values pass. |
| PD06 | complete | NLPA supplied-perpendicular, NLGC-E and NLGCE-N logarithmic solvers pass their independent A–G/selected fixture gates with residuals below 1e-8. |
| PD07 | complete | Lorentzian constant/linear and Gaussian slab kernels, exact zero-width QLT branch, both resonances, and nested refinement controls implemented. Both finite-width kernel fixtures pass. |
| PD08 | complete | Named turbulence-moment conversion, restricted balanced wave adapter and multidimensional table evaluator implemented and tested. Production provider choices remain required inputs, not defaults. |
| PD09 | partial | Analytic rigidity derivatives are available for explicit models, exact/inertial smooth QLT and NLGCE-F. Spatial gradients are supplied for constants and, when all provider gradients exist, explicit power/broken/Bohm and NLGCE-F. Integral-closure gradients, table-knot derivative policies and complete tensor divergence remain absent pending D13/host needs; absence is typed and documented. |
| PD10 | complete | Equal-length batch behavior, empty/mismatch/per-point status rules, scalar equivalence, deterministic concurrent direct calls, and changed-state sensitivity are tested. No cache is used, so stale-state reuse is impossible. |
| PD11 | deferred | D14 and host decisions remain deferred by user direction; no srcSEP3D/srcSEP binding. |
| PD12 | partial | Standalone build/tests/docs/data are packaged. Final release closure depends on PD11 and supported host build/run evidence. |

The “partial” PD09 entry does not invent derivative inputs that D13 has not
requested. A consumer needing an absent derivative must retain its coherent
provider stencil or fail the derivative-dependent operation.

## Supported inventory

All 16 first-release stable IDs are selectable:

- explicit eigenvalue models: `constant_lambda`, `constant_kappa`,
  `power_law_lambda`, `broken_rigidity_kappa`, `bohm`;
- pitch-angle/spectral models: `prescribed_lambda_mu_shape`,
  `qlt_slab_spectrum`, `qlt_slab_inertial`, `broadened_slab`;
- nonlinear pair models: `nlpa_given_perp`, `nlgc_e`, `nlgce_n`,
  `nlgce_f_2014`;
- provider/data adapters: `turbulence_adapter`, `wave_spectrum_adapter`,
  `tabulated_parallel`.

Complete SOQLT, complete composite WNLT, arbitrary directional/dynamical wave
scattering, D_pp, and non-axisymmetric perpendicular dynamics are explicitly
outside the first-release specification and are not presented as implemented.

## Acceptance evidence

| Command | Outcome |
| --- | --- |
| `cd src/models/sep_common/parallel_diffusion && make -f makefile verify` | PASS: 24 base checks plus 32 fixture-driven advanced checks. Reports: `build/test-report.json`, `build/advanced-test-report.json`. |
| `sha256sum -c parallel_diffusion_model_data/SHA256SUMS` | PASS for every bundle file. |
| `python3 parallel_diffusion_model_data/reference_verification.py` | PASS, 85 checks, 0 failed. |
| `python3 parallel_diffusion_model_data/reference_verification.py --audit` | PASS, 90 checks, 0 failed, including the 300-state audit. |
| `python3 parallel_diffusion_model_data/reference_verification.py --broadened` | PASS, 101 checks, 0 failed. |
| Strict `-Wall -Wextra -Wpedantic -Werror` `make -f makefile verify` | PASS, same 24+32 checks. |
| Address/undefined sanitizer base suite plus advanced selfcheck, NLGCE-N B, and Gaussian 0.1 fixture (`detect_leaks=0`) | PASS; no sanitizer diagnostic. |
| `cd src/models/sep_common && python3 -m json.tool SOURCE_MANIFEST.json >/dev/null && make verify` | PASS, `SEP_COMMON01`; established seven-member archive remains PIC/MPI-free with exact membership. |

The recorded bundle-verifier counts above are from the supplied verifier: the
`--audit` mode adds the audit work while retaining its 90 reported aggregate
checks; `--broadened` reports 101. Re-run commands after any numerical edit.

Backend acceptance tolerances come from `benchmark_points.json`: QLT exact
fixtures use 1e-10 relative, polynomial values 1e-11 relative, nonlinear
values 5e-8 relative with residual below 1e-8, and broadened values 1e-7
relative. The C++ advanced report records each measured discrepancy.

The advanced API selfcheck also constructs Equation (35) independently at an
`s=1` dissipation-range resonance and compares D_mu_mu at 2e-14 relative
tolerance. It separately verifies that `s_d=2` produces
`InfiniteMeanFreePath`; this is intentionally part of the aggregate advanced
selfcheck rather than an invented external fixture. The same selfcheck maps
unequal component/sign samples, cycles-per-metre data, frozen-flow frequency
data, and both Qin–Zhang reduced conventions to one canonical power, and
rejects frozen-flow input without its assumption identity.

## Remaining work

1. Resolve D13–D16 for each actual consumer, including requested derivative
   set, restart metadata, and mutual exclusion with focused D_mu_mu scattering.
2. Resolve D14 and then implement the versioned srcSEP3D `[parallel_diffusion]`
   syntax and legacy-policy decision. The current user instruction explicitly
   postpones this item.
3. Bind the shared evaluator to the existing coherent srcSEP3D coefficient
   sampling path, retain its neighbour stencil for unavailable gradients, and
   qualify scalar/batch/operator behavior in one- and four-rank native runs.
4. Only after that evidence, promote the objects into canonical `sep_common.a`
   and update both applications' exact-member/source-package audits.

No current application input selects this library, and no current result is a
claim of srcSEP3D Parker transport qualification.
