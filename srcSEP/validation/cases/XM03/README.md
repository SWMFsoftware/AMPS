# XM03 — 2013 April 11 Earth observations in Liu et al. Figure 12

XM03 validates the linked `srcSEP/AMPS` application against the Earth proton
spectra in Figure 12(a–c) of Liu et al. (2025), DOI
`10.3847/1538-4357/adc4e3`. It no longer uses SOFIE/M-FLAMPA output as the
reference and it never reads a precomputed `model/srcsep_output.csv`. Every
model spectrum is produced by the executable named with `--amps` during the
current run.

## What the case compares

`reference/liu_figure12_earth_observations.csv` contains 80 vector-extracted
measurements at 4, 12, and 36 hours after the 07:24 UTC CME flux-rope launch:

- ACE/EPAM, GOES-13/EPEAD, and SOHO/ERNE Earth measurements are included;
- STEREO-B/HET is excluded because a separate spacecraft requires a separate
  connected field line; and
- the plotted SOFIE curves are excluded because the validation target is the
  observations.

The energy-bin limits are retained as `energy_low_mev` and `energy_high_mev`.
The plotted horizontal bars describe channel width, not measurement
uncertainty. Differential intensity is in `pfu MeV^-1`, equivalent to
`(cm^2 s sr MeV)^-1`.

Liu et al. note three limitations that matter when interpreting a result: the
observations were background-subtracted, SOHO can saturate at high flux, and
ACE/EPAM does not discriminate ions. Because low-energy contamination is most
problematic, the formal score uses only points at or above 1 MeV. All points
remain visible in the CSV and comparison plots.

## Paper-derived model input

The linked calculation uses values explicitly reported by the paper:

| Quantity | XM03 value | Publication basis |
| --- | ---: | --- |
| Source-table time origin | 2013-04-11 06:00 UTC | Figure 12(d) civil-time axis |
| Model launch | 07:24 UTC (5040 s after origin) | first LASCO/C2 time and event setup |
| Earth connection | 15 min after launch | Section 4.4.2 |
| Inner/source radius | 2.5 solar radii | M-FLAMPA seed sphere |
| Earth radius | 1 AU | Table 2 |
| Earth solar wind | 363 km/s | Table 2 |
| CME/shock propagation speed | 675 km/s | EEGGL event input, Table 1 |
| Mean free path | `0.3 AU (r/AU) (pc/GeV)^(1/3)` | Equations 14–16, preferred event run |
| Diffusion | `kappa_parallel=lambda_parallel v/3` | Equation 14 |
| Injection threshold | 10 keV protons | Section 2.3.5 |
| Injection distribution | `f(p) proportional to p^-5` | Section 2.3.5 |
| Flux factor | 1.2 | postprocessing factor stated by paper |
| Perpendicular diffusion | disabled | published setup |

`input/earth_shock_thermal_source.csv` is the Earth curve from Figure 12(d).
The paper states that the injected particle number is proportional to the
plasma thermal energy density at the shock. The native callback therefore
integrates the vector trace as a piecewise-linear function and draws release
times by exact inverse-CDF sampling. It does not fit an exponential source, and
the Figure 12(a–c) observations never enter this source calculation.

## Selected heliospheric parameters and assumptions

The paper does not publish the evolving Earth-connected AWSoM field line or a
complete SOFIE restart. XM03 consequently uses the smallest transparent
one-field-line heliosphere that can run independently:

- constant 363 km/s radial wind, chosen from the event-specific Earth value in
  Table 2 rather than a generic 400 km/s solar wind;
- an equatorial Parker spiral using the Carrington 25.38-day sidereal rotation
  period (`2.8653290846e-6 rad/s`);
- the IAU nominal solar radius `6.957e8 m` and exact astronomical unit
  `149597870700 m` for unit conversion;
- a radial shock moving from 2.5 solar radii at the paper's 675 km/s EEGGL
  input speed; this is explicitly a reduced trajectory, not a reconstruction
  of the unpublished three-dimensional shock; and
- spherical expansion `div(U)=2U/r`, adiabatic momentum loss, an absorbing
  inner boundary, and first passage through 1 AU as the outward-flux sample.

The two published clocks are kept distinct: source-table hours are measured
from 06:00 UTC, while spectral snapshot hours are measured from the 07:24 UTC
flux-rope launch. The native callback adds the 5040 s offset when it samples
the 4, 12, and 36 h spectra.

The native callback uses `SEP::Transport::AdvanceParker`; the test is therefore
an application-level exercise of the production Parker SDE, spatial-diffusion
provider contract, stochastic stream, and adiabatic update. A 120 s timestep
matches the paper's reported MHD/particle coupling cadence. The 36 logarithmic
energy bins span 10 keV to 200 MeV and use 2500 particles per injection bin.
The random seed is fixed for reproducibility.

These choices do not turn the case into an exact SOFIE reproduction. The
missing global field, shock geometry, source area, and complete run deck are
listed in `publication_input.json` and must remain visible in scientific
reporting.

## Normalization and acceptance

A one-dimensional source cannot infer absolute pfu because the paper does not
publish the shock surface area or the connected flux-tube collection area.
XM03 therefore estimates exactly one multiplicative amplitude in logarithmic
space using all scored times, energies, and instruments together. It never
fits a separate factor by time, energy, or instrument. The model's spectral
and time dependence are consequently unchanged.

The provisional gates are at least 90% scored-point coverage, at most 1.0 dex
global log-RMSE, at most 0.8 dex median absolute log error, and at least 0.45
correlation between modeled and observed log intensity. The deliberately broad
thresholds recognize that this is a one-field-line Parker reconstruction, not
the paper's global MHD calculation. Tighten them only after a convergence and
observational-uncertainty study; never tune them merely to obtain PASS.

## Run

From the `srcSEP` directory:

```sh
python3 test/run_tests.py --amps ../amps \
  --validation-case XM03 \
  --output-dir test_output/XM03
```

`--case-input` is rejected for XM03 because there is one reviewed event input.
The output contains the resolved input, linked native manifest/report/log,
relative model spectrum, copied observations, point-by-point comparison CSV,
run provenance, and PNG/EPS overlays. A missing external model CSV can no
longer cause XM03 to SKIP.

Every PNG/EPS title includes the short bibliographic reference
`Liu et al. (2025), doi:10.3847/1538-4357/adc4e3` and explicitly identifies
the extracted observation panels as Figure 12(a), Figure 12(b), and Figure
12(c). The values come from `reference.plot_citation` and
`reference.figures` in `input.json`, making missing attribution a configuration
error rather than silently producing an unlabeled plot. Complete publication,
PDF/vector-figure hashes, extraction, and exclusion details remain in
`publication_input.json` and `reference/provenance.json`.

## Reproduce the vector extraction

Obtain `Fig/1304_Fig12_Spectrum_V8.pdf` from the arXiv source distribution for
`2412.07581v2`, then run:

```sh
python3 validation/reference/digitize_liu_figure12.py \
  --figure-pdf /path/to/Fig/1304_Fig12_Spectrum_V8.pdf \
  --observations validation/cases/XM03/reference/liu_figure12_earth_observations.csv \
  --earth-source validation/cases/XM03/input/earth_shock_thermal_source.csv
```

The script checks the reviewed 80-observation count. Coordinate calibration,
file hashes, exclusions, units, and limitations are recorded in
`reference/provenance.json`.
