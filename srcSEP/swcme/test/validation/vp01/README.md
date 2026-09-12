# VP01: radial solar-wind density and Leblanc profile

VP01 is a reproducible observational comparison between SWCME's normalized
Leblanc radial-density profile and the Helios 1/2 proton-core parameter archive.
It is separate from the registered synthetic `V1` regression: `V1` verifies the
metric plumbing on a small fixed fixture, while VP01 uses independently archived
spacecraft observations and produces auditable science evidence.

## Scientific question

The SWCME ambient number density is

```text
n(R) = S [3.3e5 R^-2 + 4.1e6 R^-4 + 8.0e7 R^-6] cm^-3,
```

where `R` is heliocentric radius in solar radii and one scale factor `S` is
selected so `n(1 AU)=n1AU_cm3`. VP01 asks whether the radial shape and the
frozen default normalization `n1AU_cm3=5` agree with inner-heliosphere density
observations. The outer observations estimate an equivalent 1-AU anchor for
the normalization-bias metric; that estimate is not fed back into SWCME.

The comparison does **not** assert equality of species populations. Leblanc,
Dulk, and Bougeret (1998) describe electron density, whereas the selected Helios
product reports the fitted proton core and omits the proton beam. The Helios
product documentation notes that proton-core density is typically about 10%
below the older total-proton moment product and can differ by up to 40% at
individual times. VP01 therefore validates radial shape and reports
normalization bias, not an exact electron/proton identity.

## Observational data and provenance

The data are the open **Helios corefit plasma dataset**, archived as Zenodo
record [10.5281/zenodo.1009506](https://doi.org/10.5281/zenodo.1009506) under
CC-BY-SA-4.0. Its instrument processing and scientific use are described by:

> Stansby, D., Salem, C., Matteini, L., and Horbury, T. S. (2018), “A New Inner
> Heliosphere Proton Parameter Dataset from the Helios Mission,” *Solar Physics*,
> 293, 155, [doi:10.1007/s11207-018-1377-3](https://doi.org/10.1007/s11207-018-1377-3).

The exact URLs, byte counts, MD5 values published by Zenodo, locally verified
SHA-256 values, column meanings, and species caveat are frozen in
`data/PROVENANCE.json`. Raw archives live in `data/raw/` and are ignored by Git.
`download_data.py` uses an adjacent temporary file, atomically installs a
complete download, and requires exact byte-count, MD5, and SHA-256 agreement.

Download and verify both the observations and author-supplied processing code:

```sh
cd swcme/test
make vp01-data
```

The archived `corefit.gz` file is 131,088,147 bytes (about 125 MiB compressed)
and expands logically to daily CSV files. The runner streams those files from
the tar archive; it does not extract them to disk.

## Analysis method

1. Read Helios 1 and 2 `*_corefit.csv` members directly from the verified
   archive.
2. Retain only documented status-1 fits with finite positive `n_p` and
   `0.29 <= r_sun <= 1.01 AU`.
3. Reduce every spacecraft-day to its median radius and proton-core density.
   Days, rather than individual 40-second records, receive equal weight; this
   limits artificial precision from serial correlation and cadence differences.
4. Form 14 equal-width radial bins and calculate the median and 16th/84th
   percentiles of the daily medians. A bin must contain at least 100 days.
5. Estimate an observational equivalent `n(1 AU)` by the median logarithmic
   amplitude in outer bins at `r >= 0.85 AU`, but keep the frozen SWCME default
   of 5 cm^-3 unchanged.
6. Treat every inner bin (`r < 0.85 AU`) as held-out validation of radial shape,
   and calculate the model/observed radial-slope difference over the common
   range.
7. Evaluate the bin radii twice: through the public `swcme1d::Model` C++ API and
   through `reference_solution.py`, which independently repeats the published
   coefficients and IAU distance constants.
8. Convert number density to a proton-only mass-density proxy and use SWCME's
   production thermodynamic closure to report the model mass-density change
   between zero and 5% alpha abundance. This is a composition uncertainty, not
   an alpha measurement from corefit.
9. Apply the frozen thresholds below and generate portable evidence hashes.

No raw-data value or SWCME output is manually copied into the reference
solution. The production/reference comparison can therefore detect unit,
normalization, or radial-power mistakes separately from model/observation skill.

## Acceptance criteria

| Check | Required result | Purpose |
| --- | ---: | --- |
| Accepted status-1 measurements | at least 1,000,000 | Prevent accidental partial-archive validation. |
| Daily medians | at least 2,500 | Preserve mission-scale coverage. |
| Populated radial bins | at least 12 | Require broad radial sampling. |
| C++/independent-reference maximum relative error | at most `5e-13` | Verify production implementation of the equation. |
| Median multiplicative error | at most 2.0 | Primary model/observation gate from the validation plan. |
| Absolute model/observed radial-slope difference | at most 0.30 | Primary radial-gradient gate from the validation plan. |
| Median absolute `log10(model/observation)` | at most 0.176 | Secondary factor-1.5 target. |
| 1-AU normalization bias | at most 25% | Test the frozen `n1AU_cm3` rather than fitting it away. |
| Bins inside a factor-of-two envelope | at least 80% | Require broad rather than isolated agreement. |
| Held-out median symmetric factor | at most 1.30 | Require typical inner-bin agreement within 30%. |
| Held-out worst symmetric factor | at most 1.50 | Reject a localized factor-of-two profile failure. |
| Held-out bins inside observed 16th–84th percentile | at least 80% | Compare model values with observed variability. |
| Observed binned power-law exponent | `-2.50` to `-1.80` | Require consistency with an expanding-wind trend. |

The symmetric discrepancy factor is
`max(n_model/n_observed, n_observed/n_model)`, so over- and under-prediction are
treated identically. Percentile ranges describe physical and solar-cycle
variability of daily medians; they are not independent Gaussian error bars.

## Scope and remaining release evidence

This implementation is the **Helios radial component** of VP01. A PASS means
that the checksummed Helios comparison satisfies the frozen numerical and
observational thresholds. It does not yet complete the validation plan's full
multi-mission release package. Catalog-based exclusion of ICMEs, sheaths,
shocks, and stream-interaction regions; paired PSP, Solar Orbiter, and OMNI
quiet windows; and co-located AWSoM-R, ENLIL, and EUHFORIA samples remain to be
added. `vp01_result.json` therefore records both the component PASS and
`full_multimission_release_status=INCOMPLETE` so downstream automation cannot
silently turn this first comparison into a broader science claim.

The Helios range also begins at 0.29 AU. It strongly tests the solar-wind
`R^-2` regime but has little leverage on the `R^-4` and `R^-6` coronal terms.
Those terms require radio-derived coronal density or modern inner-heliosphere
measurements and must be assessed separately before claiming full-domain
Leblanc validation.

## Running VP01

From `swcme/test`:

```sh
# Fast focused oracle test and strict C++ driver build
make vp01-test

# Run with data already present
make vp01-validation

# Or perform an end-to-end download and run
python3 validation/vp01/run_vp01.py --download
```

VP01 can also run through the global observational-validation layer. Its
`case.json` declares the shell-free launch command and result filename; all
science logic, thresholds, and data handling remain owned by this directory.

```sh
# Run just VP01 into an isolated global evidence directory.
make validation-case CASE=VP01

# Run every currently implemented observational case (VP01 today).
make validation-implemented
```

See [`../README.md`](../README.md) for the campaign registry, common artifacts,
status/exit-code semantics, and the contract used by future VP02-VP16 packages.

The Python runner requires Python 3.10 or later, NumPy, and Matplotlib. It accepts explicit
radius, binning, compiler, archive, and output-directory options; run
`python3 validation/vp01/run_vp01.py --help` for the full interface. Nonzero
exit code 1 means a scientific threshold failed. Exit code 2 means setup,
checksum, compilation, schema, or analysis failure.

## Reference and runner responsibilities

| File | Responsibility |
| --- | --- |
| `download_data.py` | Atomic retrieval and immutable checksum validation. |
| `reference_solution.py` | Independent Leblanc equation, robust outer normalization, and observed power-law fit. |
| `vp01_model_driver.cpp` | Strict-warning executable calling the public SWCME 1-D API. |
| `run_vp01.py` | Data reduction, partitioning, model/reference execution, metrics, plots, and evidence manifest. |
| `test_reference_solution.py` | Focused exact-normalization, amplitude-recovery, and slope tests. |

## Generated evidence

The default `output/` directory contains:

- `vp01_daily_medians.csv`: canonical per-day reduction;
- `vp01_binned_observations.csv`: radial medians and percentile envelopes;
- `vp01_reference_solution.csv`: independent equation results;
- `vp01_model_output.csv`: values returned by the public C++ API;
- `vp01_comparison.csv`: bin-by-bin partitions, residuals, and coverage;
- `vp01_result.json`: methods, metrics, thresholds, checks, and final status;
- `vp01_report.md`: concise human-readable result;
- `vp01_artifact_manifest.json`: byte sizes and SHA-256 hashes of evidence;
- `vp01_density_profile.png` and `.eps`: absolute profile comparison; and
- `vp01_density_ratio.png` and `.eps`: model/observation residual comparison.

The compiled driver is retained for diagnosis but excluded from the artifact
manifest because its bytes depend on compiler and platform. All portable
scientific products are hashed.

## Current frozen result

Using the checksummed Zenodo record, VP01 accepts 1,829,340 of 2,173,243 rows
and creates 3,012 daily medians in 14 bins. The outer bins imply
`n(1 AU)=4.6177 cm^-3`; the unchanged SWCME value of 5 cm^-3 has 8.28%
normalization bias. Observed and SWCME binned exponents are -2.1737 and -2.0019,
for a slope error of 0.1718. The median multiplicative error is 1.0858 and the
median absolute percentage error is 8.30%. Across the 11 held-out inner bins,
the median/worst symmetric factors are 1.0884/1.1494, and every value lies
inside the observed 16th–84th-percentile envelope. Production and independent
reference values agree to `6.78e-16`. The Helios component is **PASS**; full
multi-mission VP01 release status remains **INCOMPLETE**.
