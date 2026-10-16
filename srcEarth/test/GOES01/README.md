# GOES01: GOES-16/SGPS high-energy near-Earth transport validation

GOES01 is an end-to-end observational validation of an AMPS geospace calculation.
It launches the production Mode3D kinetic tracer, transports an independently
measured high-energy proton spectrum through a time-dependent empirical magnetic
field to the observed GOES-16 position, and compares the saved AMPS differential
intensity with GOES-16/SGPS observations. It is not a reconstruction of an SGPS
data product and it does not use GOES flux to prescribe the model boundary.

The supported result is deliberately narrow:

> For the predeclared, quality-qualified and two-head-consistent population during
> 11--13 September 2017, test whether AMPS reproduces the absolute 83--100 MeV
> proton differential intensity at geostationary orbit when driven by the observed
> SOHO/ERNE upstream spectrum and the archived T05/Qin-Denton field drivers.

This is a meaningful validation of near-Earth kinetic access/transport, angular
quadrature, unit normalization, spectral folding, time/position association, and
native output accounting. It is not a validation of SEP acceleration at the Sun or
interplanetary propagation from the Sun to 1 AU.

## Spacecraft, instrument, event, and immutable files

GOES-16 is a geostationary NOAA spacecraft located near 89.5 degrees west during
this event. The Solar and Galactic Proton Sensor (SGPS), part of the Space
Environment In-Situ Suite (SEISS), has two oppositely mounted sensor units. In the
archived no-yaw-flip geometry used here, the `-X` unit is described as west-looking
and the `+X` unit as east-looking. GOES01 names them `WEST_MINUS_X` and
`EAST_PLUS_X`; it does not infer more detailed look vectors.

The event is the 10 September 2017 solar energetic particle event associated with
the X8.2 flare/GLE72. The comparison population begins at 2017-09-11 07:30 UTC,
the interval used for the high-energy T3 comparison discussed by Kress et al.
(2021), *Observations from NOAA's Newest Solar Proton Sensor*, and ends before
2017-09-14 00:00 UTC. The first retained epoch can be later
because all observation-only qualification rules and the selected profile cadence
must also pass.

Observation root:

```text
/data/vtenishe/GEO/observational_validation_datasets_20261013/
  goes16_sgps_sep2017/
```

GOES01 authenticates these exact inputs with hashes frozen in
`goes01_common.py`:

- `se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc` is the
  NOAA special-event, five-minute, Level-2 SGPS product, version 2.0.0. The product
  contains differential proton flux, delivered uncertainty, valid-sample counts,
  data-quality sums, yaw state, and channel bounds for both sensor units.
- `spacecraft/goes16_ephemeris_1min_20170910_20170913.csv` is the derived,
  provenance-documented one-minute GOES-16 ephemeris. GOES01 uses geocentric
  latitude, east longitude, and spherical altitude from an exact position epoch;
  it never substitutes the nominal geostationary longitude.
- The data directory also preserves the original annual ephemeris NetCDF, source
  manifest, NOAA maturity/readme documents, the GOES-R product-user guide, XRS
  attitude records, and a local copy of Kress et al. (2021). Those files document
  provenance and geometry but are not silently turned into model inputs.

The corrected SGPS special-event product was generated to address documented
temperature dependence and channel-specific issues. Its flux and uncertainty units
are `protons/(cm^2 sr keV s)`. GOES01 requires those units verbatim and converts to
the AMPS per-MeV basis by exactly multiplying by 1000. It reads P8A, index 8, whose
delivered bounds differ slightly by head:

| SGPS head | P8A lower energy | P8A upper energy |
|---|---:|---:|
| `WEST_MINUS_X` (`-X`) | 83.7 MeV | 98.5 MeV |
| `EAST_PLUS_X` (`+X`) | 82.9 MeV | 99.8 MeV |

These exact bounds are verified from the NetCDF and retained in every comparison
row. AMPS is folded separately over the two intervals; the code never replaces them
with a single nominal 90 MeV point.

Model-data root:

```text
/data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926/
```

The runner accepts either that package root or its `model_inputs/SA04` directory.
It authenticates and uses:

- `validation_data/processed/SA04/boundary_soho_erne_5min.dat`: SOHO/ERNE-HED
  Level-2 proton intensity at ten nodes from 14.5 to 115 MeV, averaged to five
  minutes and converted to `m^-2 s^-1 sr^-1 MeV^-1`. It is an independent upstream
  observation. No GOES value is used to renormalize it. Package processing groups
  raw one-minute records with a floor-to-five-minute operation, and raw CDF metadata
  defines the epoch as interval start. Preparation therefore writes an auditable
  derived copy with every timestamp shifted by exactly +150 seconds and every
  intensity token unchanged. This places each average at the same interval midpoint
  used for SGPS and prevents accidental mixing of adjacent averages.
- `validation_data/processed/SA04/ts05_driver_5min.txt`: five-minute Qin-Denton
  solar-wind/geomagnetic inputs for the T05 field model.

At approximately 90 MeV a proton traverses the L1-to-near-Earth distance in much
less than one five-minute observation bin. GOES01 therefore compares the same
five-minute interval rather than applying a solar-wind convection delay, which
would be physically inappropriate for these energetic particles. It does not claim
sub-bin timing accuracy.

## Frozen observation population

Selection occurs before AMPS runs and cannot inspect a modeled value. For every
five-minute SGPS record in the fixed event interval, preparation requires:

1. exactly 300 contributing Level-1b records;
2. exactly 300 valid P8A samples in each sensor unit;
3. zero `dtc`, out-of-bounds, and error quality sums in both units;
4. finite, positive flux and supplied uncertainty in both units;
5. `YawFlipFlag == 0`;
6. an exact, good, flag-zero GOES ephemeris row at the 2.5-minute midpoint;
7. the exact midpoint-labeled ERNE average for that SGPS bin, plus real T05-driver
   rows bracketing the midpoint within the model's 900-second maximum gap; and
8. overlap of the two SGPS delivered uncertainty intervals, equivalently
   `|Jwest-Jeast| / (sigma_west + sigma_east) <= 1`.

Rule 8 is an observation-only consistency gate. The source metadata calls each
uncertainty a quadrature sum of random and systematic flux errors, but does not give
cross-head covariance. Using the sum of the two delivered uncertainties tests
interval overlap without assuming the heads' errors are independent. It demonstrates
consistency of the two broad viewing cones; it does **not** prove a fully isotropic
4-pi pitch-angle distribution.
The gate is used because exact SGPS body-frame look vectors and angular response
functions are not present in the archived package. Epochs showing measurable
east/west anisotropy are excluded rather than compared to an invented direction.
The complete ledger, including every exclusion reason, is saved as
`GOES16/input/epoch_selection_ledger.csv`.

The SGPS timestamp marks the start of a five-minute average. AMPS is evaluated at
the midpoint, and that same midpoint selects the ephemeris. This is an explicit
temporal matching rule; record-start positions are not used.

## Physical calculation

The generated deck runs `CALC_TARGET DENSITY_SPECTRUM` with the production
`GRIDDED` Mode3D solver and the HC4 mover. Despite the AMPS calculation-section
name, GOES01 validates only the saved local differential spectrum. It makes no
density comparison.

For each selected midpoint, AMPS constructs a T05 field grid and backtraces protons
from the observed GOES-16 position. The outer-boundary distribution is the
SOHO/ERNE table with logarithmic positive-intensity interpolation in time and
energy. Selected epochs coincide exactly with relabeled ERNE midpoint rows, so no
adjacent five-minute averages are blended for a comparison. The boundary spatial
distribution and pitch-angle distribution are both
explicitly uniform/isotropic. The local spectrum uses 82.9--99.8 MeV, enclosing
both P8A intervals. The comparison integrates the native
`J_local_perMeV` nodes over each head's own interval with positive power-law
interpolation and a zero-preserving linear rule at physical zeros.

The trajectory-closure policy is the same as the current SA02--SA04 production
policy. A 300-second trace is primary. Only unresolved seeds can be retried at
factor-two limits through 4800 seconds. A timeout is never relabeled trapped. The
C19-qualified bounce/drift recurrence detector must positively classify a trapped
trajectory. Any final unresolved fraction above 0.01 aborts the physical run and
cannot be converted to a nominal zero.

Mode3D uses `SNAPSHOT_LIST` and `EPOCH_PARALLEL`: each field snapshot is paired
with its one GOES location, and ranks own independent epoch meshes. OpenMP workers
parallelize the density/spectrum trajectory work inside each rank. This avoids the
incorrect epoch-by-location Cartesian product of an ordinary time series.

## Profiles

Both profiles execute independently and produce their own complete evidence and
PNG/EPS figures. They differ only in predeclared sampling resolution:

| Profile | SGPS record-start cadence | Spectrum intervals | Directions/energy | Intended use |
|---|---:|---:|---:|---|
| `QUICK` | 60 min | 24 | 96 | repository validation (`srcEarth/test/list`) |
| `COMPLETE` | 15 min | 48 | 384 | study campaign (`srcEarth/studies/list`) |

With the frozen files currently authenticated, preparation selects 54 QUICK epochs
(108 per-head comparison rows) and 217 COMPLETE epochs (434 rows). These counts are
evidence checks, not hard-coded substitutes for selection; a changed source byte is
rejected by SHA-256 before counts are computed.

## Acceptance and qualification policy

The unchanged SA01--SA04 engineering gates are applied to each head and to the
pooled per-head rows. Physical-zero AMPS values remain in every denominator.

- at least 12 matched epochs;
- maximum unresolved trajectory fraction no greater than 0.01;
- median absolute `log10(AMPS/SGPS)` no greater than 0.30 (factor two), for each
  head and pooled;
- at least 60% of each head's rows within a factor of two; and
- at least 70% of pooled rows within a factor of two.

These gates produce the **observational status**. They test a global transport
calculation against an independent boundary.
They are not SGPS calibration-accuracy limits and must not be tuned after viewing a
GOES01 residual. `comparison_summary.json` is the authoritative observational result.

A full **qualification PASS** additionally requires explicit QUICK-to-COMPLETE
numerical convergence on every per-head QUICK row. COMPLETE must contain the entire
QUICK population, use the same authenticated sources and physical/numerical contract,
and differ only in the declared refinement profile. The maximum symmetric fractional
change in the folded AMPS intensity must be no greater than 0.02. An observational
FAIL or convergence FAIL produces overall FAIL; observational PASS without linked
convergence evidence is INCOMPLETE. This prevents a visually good comparison from
concealing unresolved quadrature dependence.

`qualification_summary.json` is the authoritative overall result.
`--exit-code-policy WORKFLOW` affects only process exit status so list orchestration
can preserve plots from a physical FAIL or incomplete refinement; it never changes
the recorded status. The default `QUALIFICATION` policy returns nonzero for FAIL or
INCOMPLETE.

## Running

QUICK:

```bash
python3 srcEarth/test/GOES01/run_GOES01.py \
  --data-root /data/vtenishe/GEO/observational_validation_datasets_20261013 \
  --model-data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926/model_inputs/SA04 \
  --amps-root /data/vtenishe/GEO/AMPS --amps ./amps \
  --output-parent test_output/GOES01/quick_test_runs \
  --profile QUICK --mpi-ranks 4 --threads 6 --mover HC4 --trace-dt 0.1 \
  --formats png eps --exit-code-policy WORKFLOW
```

COMPLETE uses `--profile COMPLETE`, normally eight threads per rank, and
`test_output/GOES01/complete_study_runs`. To qualify a completed QUICK campaign, add
its exact evidence root:

```bash
python3 srcEarth/test/GOES01/run_GOES01.py ... --profile COMPLETE \
  --convergence-base test_output/GOES01/quick_test_runs/goes01_quick_<stamp>_pid<pid>
```

The runner revalidates both campaigns and writes
`comparison/quick_complete_convergence.csv` plus `convergence_summary.json` in the
COMPLETE root. A unique timestamp/PID child directory is reserved for every run;
existing physical evidence is never overwritten.

Stages `RUN`, `COMPARE`, and `PLOT` operate on an explicit `--output-root`. `ALL` and
`PREPARE` allocate a new directory below `--output-parent` unless an exact output is
provided. Preparation is fast and is useful for auditing the population without
launching AMPS:

```bash
python3 srcEarth/test/GOES01/run_GOES01.py ... --profile QUICK --stage PREPARE \
  --output-root /tmp/goes01_preparation_audit
```

No Python NetCDF package is required. `goes01_common.py` uses a small read-only
binding to the installed `libnetcdf` C runtime and verifies every requested variable
shape and required unit string. NumPy and Matplotlib are the only Python runtime
packages; see `requirements.txt`.

## Evidence and figures

Each campaign contains:

- `GOES16/prepared_run.json`: frozen selection, numerical, file-hash, and science
  contract;
- `GOES16/input/selected_observations.csv`: the two head rows used for comparison;
- `GOES16/input/epoch_selection_ledger.csv`: all event candidates and exclusions;
- `GOES16/run_status.json`, `amps.log`, and hashed native snapshot products;
- `comparison/p8a_intensity_comparison.csv` and `comparison_summary.json`;
- when linked to QUICK refinement, `comparison/quick_complete_convergence.csv` and
  `convergence_summary.json`;
- `figures/p8a_intensity_timeseries.{png,eps}`: observed per-head intensities with
  delivered uncertainties and independently transported AMPS intensities;
- `figures/p8a_model_observed_ratio.{png,eps}`: AMPS/SGPS normalization versus time;
- `figures/p8a_model_observed_scatter.{png,eps}`: model-observation scatter with
  one-to-one and factor-two lines;
- `figures/p8a_two_head_consistency.{png,eps}`: observation-only selection statistic;
- `figure_manifest.json`, `qualification_summary.json`,
  `GOES01_validation_report.md`, and
  `workflow_summary.json`.

## Limitations and prohibited interpretations

The following limitations are part of the test contract, not optional caveats:

- **No particle-density validation.** One broad high-energy differential channel
  cannot determine a proton density without assuming an unmeasured spectrum outside
  the band. GOES01 never makes that conversion.
- **No spectral-shape validation.** P8A validates one broad 83--100 MeV band. AMPS
  calculates nodes across that band, but SGPS supplies only the band product used
  here, not enough independent points to validate curvature.
- **No directional or onset validation.** Exact body-frame look vectors and angular
  response functions are absent. Head-inconsistent intervals, including strongly
  anisotropic onset behavior, are excluded by an observation-only rule.
- **Two heads are not 4-pi coverage.** Passing the z gate supports a scalar
  two-head-consistent comparison; it does not establish global isotropy.
- **Top-hat energy folding.** Only documented lower/upper channel bounds are
  available to the runner. A full energy-response matrix is not applied. Hard or
  rapidly evolving spectra can therefore introduce response-dependent bias.
- **Instrument maturity and contamination.** The special-event product contains
  corrections described in the archived NOAA and Kress documents, but SGPS has
  documented temperature, electron-contamination, background, high-energy-tail,
  and hard-spectrum caveats. The frozen quality and post-07:30 rules reduce rather
  than erase those limitations.
- **Boundary hypothesis.** SOHO/ERNE supplies the normalization and spectral shape
  at 1 AU. The boundary is modeled as spatially uniform and pitch-angle isotropic.
  GOES01 tests the consequence of that documented hypothesis; it does not validate
  it or model heliospheric SEP transport.
- **Static snapshots.** Each five-minute midpoint is a static T05 magnetic-field
  snapshot with no electric field. The test does not model field evolution during a
  particle trace or lower-energy drift/convection.
- **Field-model validity.** T05 is an empirical near-Earth model driven by the
  archived Qin-Denton series. Agreement cannot separate transport-code error from
  magnetic-field-model error; disagreement must be diagnosed, not automatically
  assigned to AMPS numerics.
- **Temporal representation.** A five-minute observed average is compared with an
  instantaneous midpoint calculation. The event selection avoids the documented
  onset interval, but sub-bin variability remains a representativeness error.
- **Single event, energy band, and orbit.** A PASS is evidence for this September
  2017 high-energy geostationary population. It does not validate lower energies,
  other geomagnetic conditions, other locations, or all AMPS geospace applications.
- **Weak-access regime.** The archived Kress study reports no evident geomagnetic
  suppression above approximately 80 MeV for its high-energy comparison interval.
  GOES01 is consequently strongest as a test of direct-access classification,
  high-energy transport bookkeeping, independent-boundary normalization, and unit/
  channel folding. It is not a stringent validation of a moving cutoff boundary;
  the existing cutoff cases retain that role.

These restrictions are why the comparison labels say “AMPS transported intensity”
and “GOES-16 SGPS observation,” never “reconstructed observation” or “perfect model.”
