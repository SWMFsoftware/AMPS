# MEAN_FREE_PATH_MODEL_DATA

Companion data for **MEAN_FREE_PATH_MODEL.md, version 1.2 (9 October 2026)**. Every transcribed number is stored as a string exactly as printed in its source, with the source key, the location in the source (table, equation, section), and a verification label. Equation numbers in these files refer to the specification.

## Layout

| Path | Content |
|---|---|
| `observations/event_fitted_mfp.csv` | Event-fitted mean free paths (Pacheco et al. 2019 Tables 2–4; Agueda et al. 2014 Table 3; Lang et al. 2024 Table 3; Lavasa et al. 2026 Tables 2–3). `quantity` is `lambda_r` (radial, SEP convention λ_r = λ∥cos²ψ) or `lambda_par`; never convert one into the other without the ψ used by the authors. |
| `parameters/sep_code_presets.json` | SEP code presets (Section 5.2). `quantity` = parallel, radial (SEP convention), isotropic scattering length, unspecified, or kappa_par. `runtime_state` and `runtime_note` give the state of Section 21.1 when the inputs named in the source are supplied. |
| `parameters/gcr_nwu_family_parameter_sets.csv` | NWU-family GCR parameter sets in long format (source, location, label, species, parameter, value_as_printed, unit, verification, flag). |
| `parameters/gcr_helmod.json` | HelMod K₀ coefficients, the numerical-use convention of D-21, and transition-function parameters. |
| `parameters/gcr_other_models.json` | Other GCR prescriptions with their normalization conventions and printed discrepancies (D-23, D-33). |
| `parameters/qlt_parameter_sets.json` | Turbulence parameter sets for the closed-form QLT expressions, the printed TS2002/TS2003 approximations, and nonlinear test-particle benchmarks. |
| `parameters/model_registry.json` | Every model ID with section, equation, source, implementation status and runtime state. |
| `turbulence/published_turbulence_values.csv` | Published turbulence quantities (slab fraction, variance, correlation lengths, spectra, model boundary values, mean field). `definition` states which of the incompatible definitions applies. |
| `benchmark_points.json` | 34 fixture groups (Section 13). `kind` says whether a value is derived, an evaluation of a printed formula, or a check of a published number. |
| `source_key_map.json` | Maps each source key to the reference identifier [Mxx] of the specification, its citation key (existing key in the supplied bibliography or proposed key in MEAN_FREE_PATH_MODEL_additions.bib), DOI and arXiv identifier; also maps equation labels to equation numbers. |
| `reference_verification.py` | Recomputes every numerical fixture value independently of the generator, compares with relative tolerances, validates the data files, runtime states and provenance fields, and checks `SHA256SUMS`. |
| `independent_physics_checks.py` | Checks the identities used by the specification against numerical quadrature, numerical differentiation and a seeded stochastic simulation. |
| `SHA256SUMS` | SHA-256 digests of all other files. |

## Verification labels

`VERIFIED-FULLTEXT xN`: read in full text, N independent readings (fetches, possibly of the same copy) agreed, counted over versions 1.0 and 1.2 of the specification; "v1.2 check" gives the readings made for version 1.2. `SECONDARY via ...`: taken from the named citing paper. `ABSTRACT-ONLY`. `FULLTEXT-GARBLED`: full text reached but layout not legible. `SIGN NOT LEGIBLE`: the sign could not be read in any accessible rendering. Flags record printed inconsistencies (Section 14 of the specification); they are never corrected silently.

## Parsing rules

A string is converted to a number only when no interpretation is needed: a plain decimal or scientific-notation number, or an exact rational literal integer/integer (for example "1/3", "5/3") in a field that is dimensionless (exponent or index). Fields containing qualifiers ("~", "≤", "≥", "+/-", ranges, "as printed"), units or unit-only normalizations are surfaced to the caller rather than converted. Normalization constants marked as units (for example "10^22 cm^2 s^-1") are not values. Derived quantities in `benchmark_points.json` are not published results.

## Checks

Requires Python 3.9+ and mpmath. Run `python3 reference_verification.py` and `python3 independent_physics_checks.py`; each exits with status 0 when every check passes. Tolerances are relative (10⁻¹² for 15-digit values; the digits printed for slopes and extremum locations); reference values that are 40-digit differences are compared with absolute bounds. The scripts check mathematics, structure and integrity; they do not check that a transcription matches its source.
